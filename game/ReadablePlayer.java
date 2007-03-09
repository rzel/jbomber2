package game;

/**
 * A <code>ReadablePlayer</code> is an object which is used to play the game,
 * but of which only data can be read, and it can not be changed
 *
 * @author <a href="mailto:m.w.manders@student.tue.nl">Maarten Manders</a>
 */
public class ReadablePlayer
{
    // The horizontal position of this Player
    private int x_pos;

    // The vertical position of this Player
    private int y_pos;

    // The id of this player
    private int id;

    // The nickname of this player
    private String nickname;

    /**
     * Creates a new <code>ReadablePlayer</code> object, and
     * assigns a given id and nickname to it.
     *
     * @param id       The id to add to this <code>ReadablePlayer</code>.
     * @param nickname The nickname of this player.
     */
    public ReadablePlayer(int id, String nickname)
    {
        this.id = id;
        this.nickname = nickname;
    }

    /**
     * Returns the horizontal position of this player.
     *
     * @return The horizontal position of this player.
     */
    public int getXPos()
    {
        return this.x_pos;
    }

    /**
     * Returns the vertical position of this player.
     *
     * @return The vertical position of this player.
     */
    public int getYPos()
    {
        return this.y_pos;
    }

    /**
     * Sets the horizontal position of this player to a new value.
     *
     * @param x_pos The new horizontal position of this player.
     */
    protected void setXPos(int x_pos)
    {
        this.x_pos = x_pos;
    }

    /**
     * Sets the vertical position of this player to a new value.
     *
     * @param y_pos The new vertical position of this player.
     */
    protected void setYPos(int y_pos)
    {
        this.y_pos = y_pos;
    }

    /**
     * Returns the id of this player.
     *
     * @return The id of this player.
     */
    public int getID()
    {
        return id;
    }

    /**
     * Returns the nickname of this player.
     *
     * @return The nickname of this player.
     */
    public String getNickname()
    {
        return nickname;
    }
}