import java.util.*;
import java.awt.event.*;

/**
 * A <code>GameEngine</code> is used by a client to handle all incoming data
 * by means of calculating player positions and other game specific things.
 * A <code>GameEngine</code> recieves data from a <Code>DelegatorClient</code>
 * and the same <code>DelegatorClient</code> is able to get data from the
 * <code>GameEngine</code>
 */
public class GameEngine
{
    // The list with players who play the same game
    private LinkedList<MovablePlayer> players;

    // The scheme that is used in this game
    Scheme scheme;

    /**
     * Creates a new <code>GameEngine</code> object.
     * The scheme that is used to play the game is specified here.
     *
     * @param scheme_name The scheme that is used to play the game.
     * @throws Exception If loading the sceme fails.
     */
    public GameEngine(String scheme_name) throws Exception
    {
        players = new LinkedList<MovablePlayer>();
        scheme = new Scheme(Config.SCHEME_PATH, scheme_name, System.out);

        // TODO initialize game engine
    }

    /**
     * Adds a new player to this game.
     *
     * id       The unique player id given by the server.
     * nickname The nickname of the player.
     * x_pos    The horizontal position of the player.
     * y_pos    The vertical position of the plater.
     */
    public void addPlayer(int id, String nickname, int x_pos, int y_pos)
    {
        // Initialize and add new player to players list
        MovablePlayer new_player = new MovablePlayer(id, nickname);
        new_player.setXPos(x_pos);
        new_player.setYPos(y_pos);
        players.addLast(new_player);
    }

    /**
     * Removes a player from the game (if it exists).
     *
     * @param id The unique player id from the player that
     *           has to be removed from the game.
     */
    public void removePlayer(int id)
    {
        // Get player with given id
        MovablePlayer mp = null;

        for (MovablePlayer p : players)
        {
            if (p.getID() == id)
            {
                // id's are unique, so only one player found
                mp = p;
            }
        }

        // If player exists, delete it from players list
        if (mp != null)
        {
            players.remove(mp);
        }
    }

    /**
     * Handles an action from a player.
     *
     * @param player_id The unique if of this player.
     * @param action_id The id of the action to handle.
     */
    public void handleAction(int player_id, int action_id)
    {
        // TODO handle moves
        MovablePlayer player = (MovablePlayer)getPlayer(player_id);

        // Player not found for some reason, return
        if (player == null) return;

        // Handle all know actions
        switch (action_id)
        {
            case Config.ACT_MOVE_RIGHT:
                player.setXPos(player.getXPos() + 1);
                break;
            case Config.ACT_MOVE_LEFT:
                player.setXPos(player.getXPos() - 1);
                break;
            case Config.ACT_MOVE_DOWN:
                player.setYPos(player.getYPos() + 1);
                break;
            case Config.ACT_MOVE_UP:
                player.setYPos(player.getYPos() - 1);
                break;
            case Config.ACT_LEAVE:
                // Remove the player that send this command
                removePlayer(player_id);
                break;
            case Config.ACT_NOTHING:
                // Do nothing
                break;
        }
    }

    /**
     * Return the player with the given id.
     *
     * @param The id of the player to find.
     * @return A <code>ReadablePlayer</code> with the given id,
     *         <code>null</code> if this player does not exist.
     */
    public ReadablePlayer getPlayer(int id)
    {
        // Search player, and return it if it exists.
        for (ReadablePlayer p : players)
        {
            if (p.getID() == id)
            {
                // id's are unique, so only one player found
                return p;
            }
        }

        return null;
    }

    /**
     * Returns an two dimensional Array that contains the scheme
     * that is used in this game.
     *
     * @return The scheme that is used in this game.
     */
    public char[][] getScheme()
    {
        return scheme.getScheme();
    }

    /**
     * Returns all the <code>ReadablePlayers</code> that play this game.
     *
     * @return all the players that play this game.
     */
    public ReadablePlayer[] getPlayers()
    {
        ReadablePlayer[] result = new ReadablePlayer[players.size()];
        // counter to keep track of players
        int i = 0;

        for (ReadablePlayer p : players)
        {
            result[i] = p;
            ++i;
        }

        return result;
    }
}